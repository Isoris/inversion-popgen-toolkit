# SPEC — Gene annotation overlay (GFF3 integration)

**Status**: drafted turn 130 final session. Not yet implemented.
Listed in scrubber roadmap (page 5 help, chat `819d8454`) as "hard" —
multi-turn work. Estimated ~2 turns total.

**Trigger** (Quentin via wishlist):
> *"Gene annotations — GFF overlay + breakpoint genes."*
>
> *"Exporter reads a GFF3 file passed via `--annotation file.gff3`,
> embeds gene coords + symbols as `genes` array in JSON. Scrubber adds
> (a) a thin gene track on page 1 below the existing tracks, (b) a
> 'genes spanned by candidate' panel on page 2, (c) a 'breakpoint-
> disrupted features' list."*

For each candidate, the atlas should show:
- Which genes overlap the inversion interval
- Which genes are disrupted by the breakpoints (5' or 3' end inside
  a gene)
- Gene density per Mb (compared to genome-wide average)

This connects inversion biology to potential **functional consequences**
— a key piece for the manuscript Discussion.

---

## 1. Three rendering surfaces

### 1.1 Page 1 — thin gene track below existing tracks

A horizontal strip showing gene boxes along the chromosome, color-coded
by orientation (+/−). Hover tooltip: `gene_symbol · biotype · length`.

Density: only show genes when zoomed in enough; collapse to a histogram
when at chromosome-wide zoom.

### 1.2 Page 2 — "Genes spanned by candidate" panel

For the focal candidate:

```
Genes overlapping inversion (n=12)
─────────────────────────────────────────────
gene_symbol  pos_start    pos_end    strand  biotype
LOC123456    15,245,123   15,267,890   +     protein_coding
SNAI2        15,512,000   15,520,123   -     protein_coding
LOC789012    15,701,200   15,709,400   +     ncRNA
...

Breakpoint-disrupted (5' end at 15.115 Mb)
─────────────────────────────────────────────
LOC456789 — disruption inside CDS (intron 4 of 7)

Breakpoint-disrupted (3' end at 17.985 Mb)
─────────────────────────────────────────────
[none — falls in intergenic region]
```

### 1.3 Page 2 — "Functional consequence" summary

A short summary block:
- Gene density: 12 / 2.87 Mb = 4.18 genes/Mb (genome-wide average:
  3.2 → above average)
- Breakpoint disruption: 1 of 2 breakpoints disrupts a gene
- Top GO terms (if InterProScan / EggNOG layer also loaded):
  *"morphogenesis (3 genes), wnt signaling (2 genes)"*

## 2. JSON layer schema

`gene_annotations_v1`:

```json
{
  "schema": "gene_annotations_v1",
  "source_gff": "fClaHyb_Gar_LG.genes.gff3",
  "genome_avg_density_per_mb": 3.2,
  "per_chromosome": {
    "C_gar_LG28": [
      {
        "gene_id": "LOC123456",
        "symbol": "SNAI2",
        "start_bp": 15245123,
        "end_bp": 15267890,
        "strand": "+",
        "biotype": "protein_coding",
        "n_exons": 7,
        "go_terms": ["GO:0009888", "GO:0007267"]
      },
      ...
    ]
  }
}
```

Atlas detects via `detectSchemaAndLayers` and computes per-candidate
gene overlap on the fly.

## 3. Producer (LANTA-side R or Python)

```r
# emit_genes_to_json.R
library(rtracklayer)
gff <- import("fClaHyb_Gar_LG.genes.gff3")
gff <- gff[gff$type == "gene"]
genes_per_chrom <- split(gff, seqnames(gff))
# Convert to JSON, write per-chromosome bundle
```

~50-100 lines.

## 4. Implementation slices

### Slice 1 — JSON layer + gene-overlap computation (~0.5 turn)
- Detect `gene_annotations_v1` layer
- `_genesOverlapping(cand)` helper
- `_genesDisruptedByBreakpoints(cand)` helper
- Tests against synthetic GFF data

### Slice 2 — page-2 gene panel (~0.5 turn)
- Render gene table per candidate
- Breakpoint-disruption summary
- Density vs genome-wide comparison

### Slice 3 — page-1 gene track (~0.5 turn)
- Thin horizontal strip canvas
- Zoom-aware rendering (genes vs density histogram)
- Hover tooltip

### Slice 4 — manuscript integration (~0.2 turn)
- Per-candidate gene list in `results_boilerplate.md`
- Supplementary table S3 (per-candidate genes spanned)

## 5. Open questions

1. **Multi-genome support**: when looking at Cgar candidates with the
   Cgar GFF vs a future Cmac analysis with the Cmac GFF — separate
   layer files per species, atlas auto-loads based on chrom prefix.
2. **GO term enrichment**: requires a reference set (genome-wide GO
   terms). Simple Fisher's exact for over-represented terms inside
   the inversion. Defer to Slice 5 if Quentin wants it.
3. **Gene biotypes to include**: protein_coding, ncRNA, tRNA, rRNA?
   All by default; user-filterable.
4. **What about pseudogenes?** Include but mark differently.

## 6. Cross-references

- `SPEC_per_candidate_breeding_readiness_card.md` — gene overlap
  feeds the breeding implication ("avoid HOM_INV × HOM_INV if one
  arrangement disrupts essential gene X").
- `SPEC_manuscript_bundle_export.md` — gene table per candidate.
- Chat `819d8454` for the original wishlist.
- Existing `SPEC_ncrna_density_layer.md` covers a sub-topic (ncRNA
  specifically).

## 7. What this is NOT

- Not gene expression analysis. We don't have RNA-seq data.
- Not mechanism inference. Whether the gene is functionally affected
  by the inversion requires functional study; the atlas only shows
  overlap.
- Not a replacement for IGV. The atlas gene track is a quick-look;
  IGV remains the right tool for read-level inspection.
