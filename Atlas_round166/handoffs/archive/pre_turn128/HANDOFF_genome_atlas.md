# Handoff — Genome Atlas (next session)

**Status**: ⚪ scaffold only · ready to start building real renderers
**File**: `Atlas/Genome_atlas.html` (27 KB, 1 page, accent orange)
**Trigger to start**: HPC under maintenance for jobs but login node is up.
The Genome Atlas's primary inputs come from the F1 hybrid assembly paper
(separate cohort, but same chromosome-scale reference). Most outputs are
already on disk from prior work.

---

## Cohort caveat — read this before doing anything

The Genome Atlas operates on **the F1 hybrid (*C. gariepinus* ×
*C. macrocephalus*) genome assembly**, NOT the 226-sample pure
*C. gariepinus* hatchery cohort. They share a chromosome-scale
reference (28 LGs) because the F1 hybrid is the assembly QC'd
to produce the reference, but they are distinct datasets with
distinct purposes:

- F1 hybrid → assembly statistics, gene tracks, repeats, synteny,
  conserved elements, ancestral karyotype reconstruction. **Genome Atlas.**
- 226 pure *C. gariepinus* → population genomics, candidate inversions,
  per-sample diversity, ROH. **Inversion Atlas + Diversity Atlas.**

The Genome Atlas is the right place to surface the assembly figures
that will appear in the F1 hybrid genome paper, plus annotations that
the Inversion Atlas needs to drill into when interpreting candidates
(gene cargo, repeat enrichment near breakpoints).

---

## Why this atlas is harder than Diversity

Three reasons:

1. **More heterogeneous data**. Diversity is essentially per-sample
   tables + per-chromosome tracks. Genome needs to handle gene
   tracks (GFF), repeat annotations (RepeatMasker .out), conserved
   elements (Cactus alignment-derived BEDs), synteny blocks
   (multi-species), and assembly stats (single big metadata blob).
2. **Cross-references to the Inversion Atlas**. Quentin's
   MODULE_CONSERVATION pipeline (VESM / SIFT4G / SnpEff / GERP++ /
   Cactus) produces deleterious-variant tracks that need to dovetail
   with candidate-zone gene cargo overlays. The data flow is bidirectional
   (Inversion Atlas → "what's at this candidate" → Genome Atlas; Genome
   Atlas → "this gene is in a known candidate" → Inversion Atlas).
3. **Manuscript figures**. Several Genome Atlas tabs map directly to
   F1 hybrid assembly paper figures. They need export modes (SVG /
   PDF / Nature Communications dimensions) that the Diversity Atlas
   doesn't.

---

## What's on disk right now (assembly side, login node accessible)

```
Genome_assembly_F1_hybrid/                     (canonical reference)
  haplotype_A.fa.gz                            ← C. gariepinus haplotype
  haplotype_B.fa.gz                            ← C. macrocephalus haplotype
  haplotype_A.gff3                             ← gene annotation A
  haplotype_B.gff3                             ← gene annotation B
  repeats/
    haplotype_A.RepeatMasker.out               ← RepeatMasker output A
    haplotype_A.repeats.bed                    ← parsed BED
    haplotype_B.RepeatMasker.out
    haplotype_B.repeats.bed
  synteny/
    Siluriformes_cactus.hal                    ← Cactus multi-species HAL
    synteny_blocks_haplotype_A_vs_*.bed        ← per-species block sets
  conserved/
    conserved_elements_haplotype_A.bed         ← Cactus-derived ≥5-way CEs
  ancestral_karyotype/
    ancestral_karyotype_blocks.bed             ← reconstructed ancestral LGs
    rearrangement_log.tsv                      ← lineage-specific events
```

Plus the WIP MODULE_CONSERVATION outputs in:

```
/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/
  module_conservation/
    vesm_650m/                                 ← VESM scores
    sift4g/                                    ← SIFT4G predictions
    snpeff/                                    ← SnpEff annotations
    gerp/                                      ← GERP++ rejected substitutions
    splice_annotation/                         ← splice-site impact
```

(Confirm exact paths next session — these are from memory and may have
been reorganized.)

---

## Recommended first cut — five tabs

### Tab 1 — Chromosome overview (highest priority)

The "ideogram bar" tab. 28 horizontal rows, one per LG, length-scaled.
Each row stacks vertically:

- centromere strip (if known) + ideogram outline
- gene density track (genes per Mb, smoothed)
- repeat density track (TE elements per Mb, optionally split by
  family — LINE / SINE / DNA / LTR)
- conserved-elements density track
- A pair of lanes for current cohort overlays:
  - L2 envelopes from active Inversion Atlas candidates (gold)
  - F_ROH heat strip from Diversity Atlas (sample × position, binned)

Click a chromosome to drill into tab 2.

This is THE manuscript figure for the F1 hybrid paper's "genome
overview" panel. Building this tab right means the manuscript figure
is essentially free — export to SVG, polish, drop in.

**Time estimate**: 6–8 hours. The R-side aggregation script is
~3 hours, the atlas renderer is ~4 hours.

### Tab 2 — Single-chromosome zoom (high priority)

Pick a chromosome, see:
- Gene track (UCSC-style, with strand glyphs)
- Repeat track (colored by class)
- Conserved-elements track
- Synteny ribbon to *Pangasius* / *Ictalurus* / outgroup
- Optional: deleterious-variant burden track from MODULE_CONSERVATION
  (if loaded)

Reads the per-chromosome slices of the same data shown in tab 1, but
at full resolution. Pannable / zoomable.

**Time estimate**: 8–10 hours. Genome browser-like UIs are not trivial.

### Tab 3 — Synteny / ancestral karyotype (manuscript-critical)

Two sub-views:

- **Block diagram**: ancestral karyotype reconstruction across
  Siluriformes, with rearrangements color-coded by lineage. Reads
  the cactus HAL summary + the reconstruction output.
- **Pairwise dot plot**: hap_A vs *C. macrocephalus* hap_B (intra-
  hybrid), or hap_A vs other Siluriformes (interspecies). Reveals
  inversions, translocations, large indels at chromosome scale.

Both views directly correspond to manuscript figures.

**Time estimate**: 10–12 hours. Dot plot rendering with sensible
binning is the hard part.

### Tab 4 — Repeat landscape (medium priority)

The TE-landscape figure: per-LG, repeat-class composition stacked-bar
plot, plus a Kimura-distance histogram showing TE age structure.
Standard genome-paper figure.

**Time estimate**: 4–5 hours.

### Tab 5 — Conserved elements + variant annotation (cross-Atlas link)

This is where MODULE_CONSERVATION outputs surface:
- Per-gene VESM / SIFT4G / SnpEff impact summary
- GERP++ conservation track
- Drill-down: click a gene, see its variant catalog

**This tab is the bridge to the Inversion Atlas**. From an Inversion
candidate's gene-cargo overlay (Inversion Atlas page 2), clicking a
gene name should jump here. Conversely, this tab needs to know which
genes fall inside any active inversion candidate (so it can
flag them visually).

**Time estimate**: 8–10 hours. Significant data volume — could be
30+ MB JSON for a full-genome variant annotation.

---

## Suggested direction — pick one

### Path A: build the chromosome overview only (tab 1) (one full day, low risk)

The single highest-value first deliverable. Produces:
- A working manuscript figure for the F1 hybrid paper
- A visual landing page that immediately shows what's in the genome
- A foundation that all other tabs share (the chromosome-row
  rendering, click-to-zoom mechanics, etc.)

The R-side script needs to:
1. Read both haplotype GFFs → emit `gene_density` per 100-kb window
2. Read both RepeatMasker outputs → emit `repeat_density`, optionally
   per-class
3. Read conserved-elements BED → emit `conserved_elements_density`
4. Optionally pull L2 envelopes from active Inversion Atlas precomps
   → emit `inversion_candidates` overlay
5. Output one JSON per haplotype: `genome_overview_haplotype_A.json`,
   `genome_overview_haplotype_B.json`

Atlas side: render the 28-row ideogram from one or both JSONs.

### Path B: build tab 5 only (deleterious variant cross-link)

Pure value for the inversion work — fills in "what genes are in this
candidate, and which variants are predicted deleterious". Skips the
manuscript-figure tabs. This is the more **useful for-Quentin's-PhD**
choice if no F1 hybrid paper push is happening this week.

Requires MODULE_CONSERVATION outputs to be in a stable shape.

### Path C: skeleton tabs 1+2+3 with empty-state styling

Don't ship real data yet — just stub each tab with a "data layer
required" panel modeled on the Inversion Atlas's page 12 / page 15
empty-state pattern. Then fill in tabs incrementally as data
sources stabilize. Lowest immediate value but lowest risk.

**My recommendation**: Path A. The chromosome overview is the highest
single-tab payoff because it doubles as a manuscript figure for the
F1 hybrid paper. After that ships, Path B (tab 5) is the natural
follow-up because it bridges to the Inversion Atlas work that's
already active.

---

## Cross-atlas navigation requirements

When the Genome Atlas is wired up to data, the following cross-atlas
clicks should work:

1. **Inversion Atlas → Genome Atlas**:
   - Click a candidate's gene-cargo entry → jump to Genome Atlas tab 5
     filtered to that gene.
   - Click a chromosome on Inversion Atlas page 1 ideogram → jump to
     Genome Atlas tab 2 zoomed to same chromosome.

2. **Genome Atlas → Inversion Atlas**:
   - On tab 1's ideogram, click a candidate-overlay strip → jump to
     Inversion Atlas page 2 with that candidate active.

3. **Genome Atlas → Diversity Atlas**:
   - On tab 1, click an F_ROH heat-strip cell → jump to Diversity Atlas
     tab 3 (chromosome ribbon) with that chromosome and sample selected.

These links are URL-based: each atlas accepts query parameters like
`?chrom=C_gar_LG07&candidate=cand_007_001` and dispatches accordingly.
Implementation note: each atlas needs a small URL-parameter-handling
function on page load.

---

## Architectural note — should Genome Atlas load Inversion Atlas
**candidate JSONs directly?**

Two options:

- **(A) Yes**: Genome Atlas reads the same per-chromosome precomp
  JSONs the Inversion Atlas reads, so tab 1 can show candidate
  overlays without any extra data plumbing. Atlas-loaded JSONs are
  cached in IndexedDB (v14), so loading once in the Inversion Atlas
  populates the cache for the Genome Atlas too — same origin.

- **(B) No**: Genome Atlas only reads genome-side JSONs. Candidate
  overlays come via URL parameter at click time, not from disk. This
  keeps the Genome Atlas conceptually genome-pure but loses the
  "always-on" overlay layer.

Recommendation: Option A. The IndexedDB cache shared across the
four atlas files makes this almost free in practice — when Quentin
loads a chromosome JSON in the Inversion Atlas, the Genome Atlas
will see it on next open. Just need to add the same
`_idbRestoreAll()` call to the Genome Atlas (and Diversity, and
Population) at page load.

---

## What this handoff is NOT

- Not a manuscript writing session — that's downstream of the figure
  rendering.
- Not for the 226 hatchery cohort's per-sample variant calls — that's
  Inversion Atlas territory.
- Not for the wild *C. macrocephalus* cohort — that's the future paper.

The Genome Atlas is for the F1 hybrid genome's chromosome-scale
assembly + annotation views, with bridges to the Inversion Atlas's
candidate-zone gene-cargo and to MODULE_CONSERVATION's variant
annotations.

---

## Open questions for next session

1. Path A, B, or C? Strong opinion: A.
2. Are the F1 hybrid GFFs and RepeatMasker outputs definitely on
   disk, or do some of those steps still need to be re-run? (The
   memory says they exist but I haven't confirmed paths.)
3. Should the Genome Atlas load Cactus HAL data directly (large,
   needs a dedicated parser) or pre-emit synteny BEDs from a
   one-shot R / python script?
4. Manuscript-figure export: SVG vs PDF? Nature Comm dimensions
   (180 mm wide for double-column)? This shapes the rendering
   layer's dimensions code.

Bring these to the next session.
