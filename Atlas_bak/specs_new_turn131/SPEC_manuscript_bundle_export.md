# SPEC — Manuscript bundle export

**Status**: drafted turn 130 final session. **Partial implementation
exists** — the page-3 catalogue already exports TSV and Markdown.
This spec documents the **bundle** version (Markdown + TSV + JSON +
figure references) that drops cleanly into a Methods/Results/Supplementary
skeleton.

**Trigger** (Quentin, multiple chats, paraphrased):
> *"Full automatic export manuscript text and send to you to write and
> finish."*
>
> *"From page 4, a 'generate bundle' button that emits a Markdown file
> (paste-ready prompt + tables + reference list) plus a parallel TSV
> (raw data) plus a JSON (the saved candidates list)."*

---

## 1. What the bundle contains

A `.zip` file with this structure:

```
manuscript_bundle_2026-05-05/
├── README.md                       # what's in here, how to use it
├── candidate_catalogue.tsv         # one row per candidate (Tier 1/2/3)
├── candidate_catalogue.json        # same data, JSON for programmatic use
├── methods_boilerplate.md          # paste-ready prose for §Methods
├── results_boilerplate.md          # paste-ready prose for §Results,
                                    # one paragraph per candidate
├── supplementary_table_S1.tsv      # full candidate list with all evidence
├── supplementary_table_S2.tsv      # per-sample karyotype matrix (226 × n_cand)
├── figures/                        # SVG exports
│   ├── F1_genome_ideogram.svg
│   ├── F2_per_candidate_grid.svg
│   └── F3_inheritance_matrix.svg
├── atlas_links.md                  # links to per-candidate atlas tabs
│                                   # for reviewers
└── generation_metadata.json        # atlas version, date, chrom set, params
```

## 2. Templated prose — what to expect

Each candidate gets a boilerplate paragraph in `results_boilerplate.md`:

> **Inversion {seq_num} — chromosome {chrom}, {start_mb}–{end_mb} Mb**
>
> A putative inversion was identified spanning {span_mb} Mb on
> chromosome {chrom}. K=3 PCA-based karyotype assignment yielded
> {n_ref} REF, {n_het} HET, and {n_inv} INV samples (n total
> = {n}). The Cochran-Armitage trend test on per-sample dosage support
> for the inverted arrangement was significant (p = {p_value};
> Bonferroni-corrected across {n_candidates} candidates).
> Sub-cluster substructure within the heterozygote class was
> {nested_or_cross_cutting} (K=6 substructure verdict: {verdict};
> {n_pure}/{n_with_data} K=6 groups pure at purity ≥ {threshold}).
>
> Inheritance Cramér's V with the next adjacent candidate
> ({adj_id} at {adj_start}–{adj_end} Mb) was {cramers_v},
> {indicating_or_not} shared chromosome ancestry. Boundary
> resolution: 5'-end {boundary_5_method} ({boundary_5_bp} ± {ci}),
> 3'-end {boundary_3_method} ({boundary_3_bp} ± {ci}).

The interpretive paragraph — what the inversion means biologically —
**is not generated**. That's Quentin's job.

## 3. What gets templated, what doesn't

### Templated (auto-fill from data)
- Genomic coordinates
- Sample counts per karyotype
- Statistical test outputs (p-values, effect sizes)
- Sub-cluster verdicts
- Cross-candidate Cramér's V
- Boundary resolution methods + uncertainties
- Cross-species synteny status (if cs_breakpoints layer present)

### NOT templated (Quentin writes)
- Biological interpretation
- Comparisons to other species' inversions
- Selection / breeding implications
- Connection to phenotype data (we don't have it)
- Discussion of the broader hatchery story

## 4. Bundle generation logic

```
function generateManuscriptBundle() {
  const candidates = state.candidateList.filter(c => c.confirmed);
  const sortedCandidates = sortByChrom(candidates);

  const bundle = new ZipBundle();

  bundle.addText('README.md', renderReadme(state, sortedCandidates));
  bundle.addText('candidate_catalogue.tsv', renderCatalogueTSV(sortedCandidates));
  bundle.addText('candidate_catalogue.json', JSON.stringify(sortedCandidates, null, 2));
  bundle.addText('methods_boilerplate.md', renderMethodsBoilerplate(state));
  bundle.addText('results_boilerplate.md',
                 sortedCandidates.map(renderCandidateParagraph).join('\n\n'));
  bundle.addText('supplementary_table_S1.tsv', renderS1(sortedCandidates));
  bundle.addText('supplementary_table_S2.tsv', renderKaryotypeMatrix(state, sortedCandidates));
  bundle.addBlob('figures/F1_genome_ideogram.svg', renderIdeogramSVG(state));
  // ... more figures
  bundle.addText('atlas_links.md', renderAtlasLinks(sortedCandidates));
  bundle.addText('generation_metadata.json', JSON.stringify({
    atlas_version: ATLAS_VERSION,
    generated_at: new Date().toISOString(),
    chrom_set: Object.keys(state.chromSummary || {}),
    params: { /* thresholds used */ },
  }, null, 2));

  return bundle.toBlob();
}
```

## 5. Implementation slices

### Slice 1 — README + catalogue TSV/JSON/MD (~0.5 turn)
- Already partially exists (page 3 catalogue export). Reuse.
- Wrap into a single `.zip` download.

### Slice 2 — boilerplate generators (~1 turn)
- `renderCandidateParagraph(cand, neighborCands)` returning Markdown string
- Methods boilerplate (one paragraph per pipeline stage that touches the
  candidate).
- Tests against synthetic candidates with known stats.

### Slice 3 — supplementary tables (~0.5 turn)
- S1: every column from candidate_catalogue.tsv plus per-evidence
  tier breakdown
- S2: 226 × n_candidate karyotype matrix

### Slice 4 — figure SVG export (~1 turn)
- Genome-wide ideogram (depends on `SPEC_genome_wide_ideogram.md`)
- Per-candidate composite (PCA + dosage strip + tier badge)
- Inheritance Cramér's V matrix

### Slice 5 — atlas link generation (~0.3 turn)
- For each candidate, generate `atlas_links.md` entries pointing reviewers
  at the relevant atlas page anchors.

## 6. Open design questions

1. **Where does the button live?** Probably page-3 catalogue toolbar,
   alongside the existing TSV/MD exports.
2. **How heavy can the bundle get?** With 100+ confirmed candidates
   plus per-candidate SVG figures, easily 10+ MB. Browser needs to
   handle that.
3. **Per-candidate folders inside `figures/`?** Probably yes for ≥10
   candidates. Subfolder per candidate with composite + per-track
   plots.
4. **Citation handling**: Quentin's bibliography lives elsewhere. The
   bundle should emit citation keys (e.g., `\cite{Andres2026}`) but not
   the actual references.
5. **Language**: English-only. Thai i18n is a separate concern (already
   spec'd in `SCHEMA §18`).

## 7. Honest expectations

This bundle replaces ~2 days of manuscript prep work per round of
candidates. It does NOT replace:

- Reading the candidates and writing the discussion
- Justifying tier thresholds in Methods
- Drawing per-candidate biological conclusions
- Writing the Abstract / Significance statement

What it provides: numbers Quentin can trust, in a form that drops into
the manuscript skeleton without manual TSV-wrangling.

## 8. Tests

- Generate bundle from 5 synthetic confirmed candidates → all 8 files
  present, valid Markdown / TSV / JSON.
- Boilerplate paragraph contains all expected fields filled in (no
  `{placeholder}` strings remain).
- Karyotype matrix has correct dimensions (226 × 5).
- Supplementary table rows match candidate count.

## 9. Dependencies

- `state.candidateList` with `confirmed` candidates (existing).
- `SPEC_genome_wide_ideogram.md` for the ideogram SVG.
- `SPEC_per_candidate_breeding_readiness_card.md` for per-candidate
  composites.
- Existing TSV/MD export logic (page 3).
