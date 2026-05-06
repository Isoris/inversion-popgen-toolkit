# SPEC — Per-candidate breeding-readiness card (Atlas 5 Part B)

**Status**: drafted turn 130 final session. Spec'd in chat `c03fc41e`
("Plan arrangement and review") as Atlas 5 Part B. Not yet implemented.

**Trigger** (Quentin, chat `c03fc41e`):
> *"Atlas 5 — Per-individual integration atlas / Per-candidate
> breeding-readiness atlas (Supplementary). [...] Part B — Per-candidate
> breeding-readiness card (one page per marker-ready candidate)."*
>
> *"Atlas 5 is the deliverable that converts the paper from 'a population
> genomics study' into 'a hatchery management resource.'"*

A printable per-candidate one-pager designed for assay developers and
hatchery geneticists who don't read the full manuscript but need to
make decisions about marker design or pairing strategy.

---

## 1. Card contents (one candidate per page)

```
╔═══════════════════════════════════════════════════════════════════╗
║  Inversion {seq_num} — chromosome {chrom}                       ║
║  Coordinates: {start_mb}–{end_mb} Mb (span {span_mb} Mb)        ║
║  Tier: {tier}   |  Confidence: {confidence}                     ║
╠═══════════════════════════════════════════════════════════════════╣
║                                                                   ║
║  [PCA panel — 226 samples colored by REF/HET/INV]                ║
║                                                                   ║
║  Karyotype counts:  REF: {n_ref}   HET: {n_het}   INV: {n_inv}   ║
║                                                                   ║
║  Carrier-by-ancestry table:                                       ║
║   K=8 cluster    | n REF | n HET | n INV | total                 ║
║   broodline 1    |  12   |  8    |  2    |   22                  ║
║   broodline 2    |  4    |  18   |  16   |   38                  ║
║   ...                                                             ║
║                                                                   ║
║  Per-arrangement burden:                                          ║
║   REF mean F_ROH: {fr_ref}   damaging-load: {dl_ref}             ║
║   INV mean F_ROH: {fr_inv}   damaging-load: {dl_inv}             ║
║   F_ROH differential (Wilcoxon): p = {p_value}                   ║
║                                                                   ║
║  Recombinant fraction: {n_recom} / {n_het} HETs ({pct}%)         ║
║                                                                   ║
║  Marker panel: see `marker_panel_inv{N}.tsv`                      ║
║   - {n_markers} private markers (Tier 1)                          ║
║   - {n_panel} multi-marker panels (≥3 linked Tier 1)              ║
║   - PCR feasibility: {n_easy} easy, {n_hard} hard                 ║
║                                                                   ║
║  Pairing-design implication:                                      ║
║   {auto_pairing_advice}                                           ║
║                                                                   ║
║  Cross-species orthology:                                         ║
║   - Cgar reference: present (this study)                          ║
║   - Cmac reference: {synteny_status}                              ║
║   - Other Clarias: {n_other_species_supporting}                   ║
║                                                                   ║
║  Atlas reviewer link: <atlas-url>#cand={candidate_id}             ║
║                                                                   ║
╚═══════════════════════════════════════════════════════════════════╝
```

## 2. Auto-pairing advice logic

Templated based on the per-arrangement burden differential:

| Condition | Auto-advice text |
|---|---|
| INV F_ROH significantly higher (p < 0.05) | "Avoid HOM_INV × HOM_INV pairings — INV-arrangement carriers show elevated F_ROH ({delta} above REF mean). Cross HOM_INV with HOM_REF or HET to break ROH." |
| Damaging-load asymmetric | "INV arrangement carries {n_extra_dam} additional predicted damaging variants per Mb relative to REF. Prefer HOM_REF as paternal contributor." |
| MAF imbalanced (one arrangement near fixation) | "Arrangement frequencies skewed ({maf}). Consider expanding broodstock from sources carrying the minor arrangement." |
| Recombinant fraction high (>20% of HETs) | "Recombinant background substantial. PCR markers should target arrangement-distinguishing sites, not boundary regions." |
| Otherwise | "No specific pairing constraint indicated; standard MAF-balanced pairing applies." |

These are **suggestions**, not hard prescriptions. The card explicitly
flags them as "auto-generated from genomic statistics; final breeding
decisions require human review."

## 3. Generation pipeline

Atlas-side function:

```
function generateBreedingCard(cand) {
  return {
    chrom: cand.chrom,
    coordinates: { start: cand.start_bp, end: cand.end_bp, span_mb: ... },
    tier: cand.tier,
    confidence: cand.confidence,
    pca_panel_svg: renderPCAPanel(cand),
    karyotype_counts: countKaryotypesByBand(cand),
    ancestry_table: crossTabKaryotypeByQ8(cand, state.q8),
    burden: computePerArrangementBurden(cand, state.fr_oh, state.damaging),
    recombinant_fraction: computeRecombinantFraction(cand),
    marker_panel: state.markerPanels[cand.id] || null,
    auto_advice: generatePairingAdvice(cand),
    cross_species: lookupSynteny(cand),
    atlas_url: `${ATLAS_BASE}#cand=${cand.id}`,
  };
}
```

Output: HTML or PDF page per candidate (use the existing
print-ready-CSS pattern from page 2 candidate focus).

## 4. Implementation slices

### Slice 1 — card data assembly (~0.5 turn)
- `generateBreedingCard(cand)` data builder
- Auto-advice text generator with the rules above
- Tests against synthetic candidates

### Slice 2 — print-ready HTML render (~0.5 turn)
- One-page HTML layout matching the visual contract above
- Inline SVG for PCA panel
- Print stylesheet for `@page` size A4

### Slice 3 — bulk PDF export (~0.5 turn)
- "Generate breeding cards" button on the catalogue page
- Per-tier filter (Tier 1/2 by default)
- Output: zip of PDFs, one per candidate

## 5. Dependencies

- Confirmed candidate list with locked karyotype labels.
- Per-sample F_ROH (from `sample_froh` JSON layer).
- Per-sample K=8 ancestry assignment (from NGSadmix Q output).
- Per-sample damaging-load summary (from `MODULE_CONSERVATION` output).
- Cross-species synteny status (from `cs_breakpoints_v1` layer).
- Marker panel TSV (from `SPEC_marker_panel_design_atlas.md`).

## 6. What this is NOT

- Not a phenotype-linked pairing optimizer. We don't have phenotype data.
- Not an automated breeding-program planner. It's a per-candidate
  diagnostic, intended to feed into breeding decisions made by human
  experts.
- Not a replacement for the breeder's own broodstock records. It only
  reflects what the genomic data shows.

## 7. Tests

- Synthetic candidate with planted F_ROH asymmetry → auto-advice
  triggers the right rule.
- Multi-candidate batch generation → all cards rendered, no missing
  fields.
- Print stylesheet → A4 layout fits content.

## 8. Cross-references

- `SPEC_marker_panel_design_atlas.md` for the marker-panel column.
- `SPEC_recombinant_dosage_changepoint_detector.md` for the recombinant
  fraction.
- `SPEC_manuscript_bundle_export.md` for inclusion in Supplementary
  Data export.
- Original spec lives in chat `c03fc41e` "Plan arrangement and review"
  Atlas 5 Part B section.
